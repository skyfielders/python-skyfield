
all:
	echo "Targets: fonts, novas, spice"

# Download and unzip:
# http://scripts.sil.org/cms/scripts/render_download.php?format=file&media_id=GentiumPlus-5.000-web.zip&filename=GentiumPlus-5.000-web.zip

GENTIUM_IN := $(HOME)/Downloads/GentiumPlus-5.000-web/web/GentiumPlus-R.woff
WOFF_OUT := ../skyfield/documentation/_static/GentiumPlus-R-custom.woff
TEXTS := $(wildcard ../skyfield/documentation/*.rst)

fonts: $(GENTIUM) $(TEXTS)
	python2.7 strip_woff.py $(GENTIUM_IN) $(WOFF_OUT) $(TEXTS)
	ls -l $(GENTIUM_IN)
	ls -l $(WOFF_OUT)

spice: ../skyfield/data/spice.py

../skyfield/data/spice.py: build_spice.py ~/spice/toolkit/src/spicelib/chgirf.f
	python build_spice.py ~/spice/toolkit > $@

novas: ../skyfield/tests/test_against_novas.py

../skyfield/tests/test_against_novas.py: build_novas_tests.py
	python $< > tmp.py && mv tmp.py $@
