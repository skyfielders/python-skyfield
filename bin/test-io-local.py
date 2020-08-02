#!/usr/bin/env python

import os
import threading
from contextlib import contextmanager

try:
    from http.server import HTTPServer, SimpleHTTPRequestHandler
except ImportError:
    from BaseHTTPServer import HTTPServer
    from SimpleHTTPServer import SimpleHTTPRequestHandler

from skyfield.api import load

@contextmanager
def cd(directory):
    old = os.getcwd()
    try:
        os.chdir(directory)
        yield
    finally:
        os.chdir(old)

@contextmanager
def local_web_server():
    httpd = HTTPServer(('localhost', 0), SimpleHTTPRequestHandler)
    thread = threading.Thread(target=httpd.serve_forever, daemon=True)
    thread.start()
    yield httpd.server_address

def main():
    here = os.path.dirname(__file__)
    with cd(os.path.join(here, '..', 'ci')):
        with local_web_server() as (host, port):
            load('http://{0}:{1}/de421.bsp'.format(host, port), reload=True)

if __name__ == '__main__':
    main()
