# -*- coding: utf-8 -*-
"""Parsing routines for JPL text PCK files.

For an (incomplete) summary of the file format, look for the heading
“NAIF Text Kernel Format” in the “PCK Required Reading”:

https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html

"""
import re

def load(lines, variables):
    """Load PCK text kernel names and values into a ``variables`` dict."""
    for name, equals, values in parse(lines):
        if equals == b'=':
            if len(values) == 1:
                variables[name] = values[0]  # unbox a scalar value
            else:
                variables[name] = values
        elif equals == b'+=':
            old = variables.get(name)
            if old is None:
                old = variables[name] = []
            elif not isinstance(old, list):
                old = [old]  # box an earlier scalar before extending
            old.extend(values)

def parse(lines):
    """Yield ``(name, equals, values)`` bytestrings from a PCK text kernel.

    This merely reads raw assignment statements; it doesn’t combine
    multiple ``+=`` assignments to create single values.  The byte
    string ``equals`` will be either ``b'='`` or ``b'+='``.  Scalars are
    returned as a ``values`` list one item long.

    """
    tokens = iter(_parse_tokens(lines))
    for token in tokens:
        name = token.decode('ascii')
        equals = next(tokens)
        if equals not in (b'=', b'+='):
            raise ValueError('an equals sign is expected after %r' % name)
        token = next(tokens)
        if token == b'(':
            values = []
            for token in tokens:
                if token == b')':
                    break
                values.append(_evaluate(token))
        else:
            values = [_evaluate(token)]
        yield name, equals, values

def _evaluate(token):
    """Return a string, integer, or float parsed from a PCK text kernel."""
    if token[0:1].startswith(b"'"):
        return token[1:-1].decode('ascii')
    if token.isdigit():
        return int(token)
    if token.startswith(b'@'):
        raise NotImplementedError('TODO: need parser for dates,'
                                  ' like @01-MAY-1991/16:25')
    token = token.replace(b'D', b'E')  # for numbers like -1.4D-12
    return float(token)

_token_re = re.compile(b"[A-Za-z]\\w+|=|\\+=|\\(|\\)|'[^']*'|[^), ]+")

def _parse_tokens(lines):
    """Yield all the tokens inside the data segments of a PCK text file."""
    lines = iter(lines)
    for line in lines:
        if b'\\begindata' not in line:
            continue  # save cost of strip() on most lines
        line = line.strip()
        if line != b'\\begindata':
            continue
        for line in lines:
            line = line.strip()
            if line == b'\\begintext':
                break
            for token in _token_re.findall(line):
                yield token
