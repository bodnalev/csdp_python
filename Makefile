all:
	python3 setup.py build_ext --inplace

clean:
	rm -rf src/*.out src/*.bin src/*.exe src/*.o *.a *.so test build