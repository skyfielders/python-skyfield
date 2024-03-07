
def setup(app):
    app.connect('autodoc-process-signature', signature)

def signature(app, what, name, obj, options, signature, return_annotation):
    #print(repr(signature))
    pass
