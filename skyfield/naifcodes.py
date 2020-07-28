from jplephem.names import (
    target_name_pairs as code_name_pairs,
    target_names as code_names
)

name_codes = dict((name, code) for code, name in code_name_pairs)

def numbered_name_of(code):
    """Given a code, return a string giving both the code and name.

    >>> numbered_name_of(301)
    '301 Moon'

    """
    name = code_names.get(code, '(Unnamed)')
    return '{0} {1}'.format(code, name)

def _target_name(target):
    """Return `target` annotated for display to the user.

    * A string target is quoted to make clear that it's a string, like 'comet'.
    * A numeric target has its NAIF name appended, like 399 EARTH.

    """
    if isinstance(target, str):
        return repr(target)
    name = code_names.get(target, 'UNKNOWN')
    return '{0} {1}'.format(target, name)
