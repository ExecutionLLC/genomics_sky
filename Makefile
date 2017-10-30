all:
		python setup.py build_ext --inplace

sdist:
		python setup.py sdist

test:		all
		python runtests.py

clean:
		@echo Cleaning Source
		@rm -rf build
		@rm -rf dist
		@rm -f *.py[cdo] */*.py[cdo] */*/*.py[cdo] */*/*/*.py[cdo]
		@rm -f *.so */*.so */*/*.so
		@rm -f *~ */*~ */*/*~
