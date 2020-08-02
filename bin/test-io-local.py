#!/usr/bin/env python

import os
import tempfile
import threading
from contextlib import contextmanager

try:
    from http.server import HTTPServer, SimpleHTTPRequestHandler
except ImportError:
    from BaseHTTPServer import HTTPServer
    from SimpleHTTPServer import SimpleHTTPRequestHandler

from skyfield.api import Loader

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
    directory = tempfile.mkdtemp()
    load = Loader(directory)
    here = os.path.dirname(__file__)
    with cd(os.path.join(here, '..', 'ci')):
        with local_web_server() as (host, port):
            # The first time, it should create the file in the empty directory.
            load('http://{0}:{1}/de421.bsp'.format(host, port), reload=True)
            filenames = os.listdir(directory)
            print(filenames)
            assert filenames == ['de421.bsp']

            # The second time, it should overwrite it.
            load('http://{0}:{1}/de421.bsp'.format(host, port), reload=True)
            filenames = os.listdir(directory)
            print(filenames)
            assert filenames == ['de421.bsp']

            # The third time, it should make a backup.
            load('http://{0}:{1}/de421.bsp'.format(host, port), reload=True,
                 backup=True)
            filenames = os.listdir(directory)
            print(filenames)
            assert filenames == ['de421.bsp', 'de421.old1.bsp']

if __name__ == '__main__':
    main()
