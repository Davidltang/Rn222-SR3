"""
Execute code on a remote Jupyter kernel via websocket.
Connects to the Midway3 JupyterLab kernel through the SSH tunnel.
"""
import sys
import json
import requests
from jupyter_client import BlockingKernelClient

# Auth setup from the cookie jar
import http.cookiejar
cookie_jar = http.cookiejar.MozillaCookieJar('/tmp/jupyter_cookies_new')
try:
    cookie_jar.load()
except:
    pass

BASE = "http://localhost:16812"
KERNEL_ID = "8a4e4ffc-0f1a-4af2-ad68-f6dddbc5d1b9"

# Get XSRF token from cookies
xsrf = None
for cookie in cookie_jar:
    if cookie.name == '_xsrf':
        xsrf = cookie.value
        break

if not xsrf:
    print("ERROR: No XSRF token found in cookies")
    sys.exit(1)

session = requests.Session()
session.cookies = cookie_jar
session.headers.update({"X-XSRFToken": xsrf})

# Get kernel info to verify
resp = session.get(f"{BASE}/api/kernels/{KERNEL_ID}")
if resp.status_code != 200:
    print(f"ERROR: Cannot access kernel: {resp.status_code}")
    sys.exit(1)

kinfo = resp.json()
print(f"Connected to kernel: {kinfo['id'][:12]}... ({kinfo['execution_state']})")

# Connect via websocket using the kernel client
ws_url = f"ws://localhost:16812/api/kernels/{KERNEL_ID}/channels"
# Add cookie and XSRF token as query params for websocket auth
ws_url_with_auth = f"{ws_url}?token={xsrf}"

print(f"WebSocket URL: {ws_url}")

# Use jupyter_client's client
kc = BlockingKernelClient()
kc.load_connection_info({
    "transport": "tcp",
    "ip": "localhost",
    "kernel_name": "python3",
    "shell_port": 0,
    "iopub_port": 0,
    "stdin_port": 0,
    "hb_port": 0,
    "control_port": 0,
})

# We need to manually set up the websocket connection
# jupyter_client doesn't directly support remote websocket URLs out of the box
# Let's use raw websocket instead

print("\nUsing raw websocket approach...")

import websocket
import threading
import uuid
import time

# Cookie string for websocket connection
cookie_str = "; ".join(f"{c.name}={c.value}" for c in cookie_jar)

class KernelExecutor:
    def __init__(self, ws_url, cookies, xsrf_token):
        self.ws_url = ws_url
        self.cookies = cookies
        self.xsrf_token = xsrf_token
        self.ws = None
        self.results = {}
        self.iopub_messages = []
        self.lock = threading.Lock()
        self.msg_id = None

    def connect(self):
        headers = {
            "Cookie": self.cookies,
        }
        self.ws = websocket.create_connection(
            self.ws_url,
            header=headers,
            timeout=30
        )
        print("WebSocket connected!")

    def _on_message(self, ws, message):
        """This is not used for sync websocket - we poll instead."""
        pass

    def execute(self, code, timeout=60):
        """Execute code and return the result."""
        import datetime

        msg_id = str(uuid.uuid4())
        session_id = str(uuid.uuid4())[:8]

        # Create execute_request message
        header = {
            "msg_id": msg_id,
            "msg_type": "execute_request",
            "username": "claude",
            "session": session_id,
            "date": datetime.datetime.now().isoformat(),
            "version": "5.3",
        }

        content = {
            "code": code,
            "silent": False,
            "store_history": True,
            "user_expressions": {},
            "allow_stdin": False,
            "stop_on_error": True,
        }

        parent_header = {}
        metadata = {}

        msg = {
            "header": header,
            "parent_header": parent_header,
            "metadata": metadata,
            "content": content,
            "channel": "shell",
            "buffers": [],
        }

        # Send
        self.ws.send(json.dumps(msg))

        # Collect response
        stdout_text = []
        stderr_text = []
        execute_result = None
        error_result = None
        status = "unknown"

        start = time.time()
        while time.time() - start < timeout:
            try:
                self.ws.settimeout(1.0)
                raw = self.ws.recv()
                resp = json.loads(raw)

                msg_type = resp.get("header", {}).get("msg_type", "")
                parent_id = resp.get("parent_header", {}).get("msg_id", "")

                # Only process messages related to our request
                if parent_id == msg_id or msg_type == "status":
                    if msg_type == "stream":
                        name = resp.get("content", {}).get("name", "stdout")
                        text = resp.get("content", {}).get("text", "")
                        if name == "stdout":
                            stdout_text.append(text)
                        elif name == "stderr":
                            stderr_text.append(text)
                    elif msg_type == "execute_result":
                        execute_result = resp.get("content", {}).get("data", {})
                    elif msg_type == "error":
                        error_result = resp
                    elif msg_type == "status":
                        status = resp.get("content", {}).get("execution_state", status)
                    elif msg_type == "execute_reply":
                        status = resp.get("content", {}).get("status", status)
                        break  # execute_reply is the final response

            except websocket.TimeoutError:
                continue
            except Exception as e:
                print(f"WebSocket error: {e}")
                break

        return {
            "status": status,
            "stdout": "".join(stdout_text),
            "stderr": "".join(stderr_text),
            "result": execute_result,
            "error": error_result,
        }

    def close(self):
        if self.ws:
            self.ws.close()

# Execute
executor = KernelExecutor(ws_url, cookie_str, xsrf)
try:
    executor.connect()

    # Simple test
    print("\n--- Checking kernel state ---")
    result = executor.execute("""
print("Hello from kernel!")
print(f"df exists: {'df' in dir()}")
if 'df' in dir():
    print(f"df shape: {df.shape}")
    print(f"df columns: {list(df.columns)}")
print(f"df_peaks exists: {'df_peaks' in dir()}")
if 'df_peaks' in dir():
    print(f"df_peaks shape: {df_peaks.shape}")
    print(f"df_peaks columns: {list(df_peaks.columns)}")
print(f"df_fv exists: {'df_fv' in dir()}")
if 'df_fv' in dir():
    print(f"df_fv shape: {df_fv.shape}")
    has_rn_band = 'in_rn_band' in df_fv.columns
    print(f"df_fv has in_rn_band col: {has_rn_band}")
print(f"p_opt_rn exists: {'p_opt_rn' in dir()}")
if 'p_opt_rn' in dir():
    print(f"p_opt_rn: {p_opt_rn}")
print(f"exp_bump_model exists: {'exp_bump_model' in dir()}")
""")

    print(f"Status: {result['status']}")
    print(f"STDOUT:\n{result['stdout']}")
    if result['stderr']:
        print(f"STDERR:\n{result['stderr']}")
    if result['error']:
        print(f"ERROR:\n{result['error']}")

finally:
    executor.close()
