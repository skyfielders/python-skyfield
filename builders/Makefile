
all:
	echo "Targets: arrays, fonts, magnitudes, novas, spice"

arrays:
	cd .. && rm -f finals2000A.all && python builders/build_arrays.py

# Download and unzip:
# http://scripts.sil.org/cms/scripts/render_download.php?format=file&media_id=GentiumPlus-5.000-web.zip&filename=GentiumPlus-5.000-web.zip

GENTIUM_OUT := ../skyfield/documentation/_static/GentiumPlus-R-custom.woff
UBUNTU_OUT := ../skyfield/documentation/_static/UbuntuMono-R-custom.woff

fonts: $(GENTIUM_OUT) $(UBUNTU_OUT)

GENTIUM_IN := $(HOME)/Downloads/GentiumPlus-5.000-web/web/GentiumPlus-R.woff
UBUNTU_IN := /usr/share/fonts/truetype/ubuntu-font-family/UbuntuMono-R.ttf

TEXTS := $(wildcard ../skyfield/documentation/*.rst)

$(GENTIUM_OUT): $(GENTIUM_IN) $(TEXTS)
	python2.7 strip_woff.py $(GENTIUM_IN) $(GENTIUM_OUT) $(TEXTS)
	ls -l $(GENTIUM_IN)
	ls -l $(GENTIUM_OUT)

$(UBUNTU_OUT): $(GENTIUM_IN)
	python2.7 strip_woff.py $(UBUNTU_IN) $(UBUNTU_OUT) ascii.txt
	ls -l $(UBUNTU_IN)
	ls -l $(UBUNTU_OUT)

spice: ../skyfield/data/spice.py

../skyfield/data/spice.py: build_spice.py ~/spice/toolkit/src/spicelib/chgirf.f
	python build_spice.py ~/spice/toolkit > $@

magnitudes: ../skyfield/tests/test_magnitudes_raw.py

../skyfield/tests/test_magnitudes_raw.py: build_magnitude_tests.py Ap_Mag_Input_V3.txt
	python build_magnitude_tests.py > $@
	cat $@

novas: ../skyfield/tests/test_against_novas.py

../skyfield/tests/test_against_novas.py: build_novas_tests.py
	python $< > tmp.py && mv tmp.py $@

DE430 := https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
JUP310 := https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup310.bsp

excerpts: \
 ../skyfield/tests/data/de430-2015-03-02.bsp \
 ../skyfield/tests/data/jup310-2015-03-02.bsp \
 ../skyfield/tests/data/jup310-2053-10-08.bsp \
 ../ci/hip_main.dat.gz

../skyfield/tests/data/de430-2015-03-02.bsp:
	python -m jplephem excerpt 2015/3/2 2015/3/3 $(DE430) $@

../skyfield/tests/data/jup310-2015-03-02.bsp:
	python -m jplephem excerpt 2015/3/2 2015/3/3 $(JUP310) $@

../skyfield/tests/data/jup310-2053-10-08.bsp:
	python -m jplephem excerpt 2053/10/8 2053/10/9 $(JUP310) $@

../ci/hip_main.dat.gz: ../hip_main.dat.gz
	zcat $< | awk -F\| '$$6 <= 6.6 || $$2 == 87937' | gzip -c -9 > $@
