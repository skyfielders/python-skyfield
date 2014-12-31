from skyfield.units import Angle

def new_cross(catalog, name):
    '''Uses the New Cross Index (`IV/27A <http://cdsarc.u-strasbg.fr/viz-bin/Cat?IV/27A>`_) to look up designations between catalogs. It also includes the visual magnitude, right ascension, magnitude and constellation abbreviation.

    Designations:

    Henry Draper Catalog Number ``HD`` (`III/135 <http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=III/135>`_),
    Durchmusterung Identification from HD Catalog ``DM`` (see `IV/27A note 1 <http://cdsarc.u-strasbg.fr/viz-bin/Cat?IV/27A#sRM3.1>`_),
    General Catalogue of 33342 stars ``GC`` (`I/113 <http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/113>`_),
    Harvard Revised Number - BSC5 ``HR`` (`V/50 <http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=V/50>`_),
    Hipparcos Catalog ``HIP`` (`I/196 <http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=I/196>`_)

    You can look solely for the Bayer ``Bayer`` *or* Flamsteed number ``Fl``, or return either with ``BFD``.

    `IV/27A note 2 <http://cdsarc.u-strasbg.fr/viz-bin/Cat?IV/27A#sRM3.2>`_ on HIP and CSI (`IV/9 <http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=IV/9>`_) right ascensions and visual magnitude:

        Right ascensions, declinations and visual magnitudes for all stars were taken from the Hipparcos catalog and from the CSI for the stars that has no number in catalog Hipparcos.

    Here is an example of looking up a star in Orion's belt by using the Hipparcos number:

    ::

        from skyfield.data import new_cross
        mintaka_names = new_cross.new_cross('HIP','25930')
        print(mintaka_names)

    Excluding DE to Dec as a key abbreviation, the labels are the same as the New Cross Index.
    '''

    cross_index_url = 'ftp://cdsarc.u-strasbg.fr/pub/cats/IV/27A/catalog.dat'
    star_dict = {}

    try:
        data = open('catalog.dat','r')
    except:
        import wget
        wget.download(cross_index_url)
        data = open('catalog.dat','r')

    for l in data.readlines():
        s = {}
        s['HD']   = l[0:6]
        s['DM']   = l[7:19]
        s['GC']   = l[20:25]
        s['HR']   = l[26:30]
        if(s['HR'] == '    '): s['HR'] = None
        s['HIP']  = l[31:37]
        if(s['HIP'] == '      '): s['HIP'] = None
        ra = float(l[38:40])+float(l[40:42])/60.+float(l[42:47])/3600.
        s['RA']   = Angle(degrees=float(ra))
        if(l[48]=='+'):
            sign=1
        else:
            sign=-1
        dec = sign*float(l[49:51])+float(l[51:53])/60.+float(l[53:57])/3600.
        s['Dec']  = Angle(degrees=float(dec))
        s['Vmag'] = l[58:63]
        s['Fl']   = l[64:67]
        s['Bayer']= l[68:73]
        if(s['Bayer'] == '     '):
            s['Bayer'] = None
        if(s['Bayer']):
            s['BFD'] = s['Bayer']
        elif(s['Fl']):
            s['BFD'] = s['Fl']
        else:
            s['BFD'] = None
        s['Cst']  = l[74:77]
        if(s[catalog] is not None and s[catalog].lower().strip() == name.lower().strip()):
            return s
    return False
