
all: ../skyfield/tests/test_against_novas.py

../skyfield/tests/test_against_novas.py: build_novas_tests.py
	python $< > tmp.py && mv tmp.py $@



horizons_output.txt: horizons_input.txt
	python horizons_telnet.py
