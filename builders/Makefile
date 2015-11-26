
all: ../skyfield/data/spice.py

../skyfield/data/spice.py: build_spice.py ~/spice/toolkit/src/spicelib/chgirf.f
	python build_spice.py ~/spice/toolkit > $@

all: ../skyfield/tests/test_against_novas.py

../skyfield/tests/test_against_novas.py: build_novas_tests.py
	python $< > tmp.py && mv tmp.py $@
