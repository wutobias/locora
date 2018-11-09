install:
	python setup.py install

clean:
	rm -rf *.egg-info
	rm -rf build
	rm -rf dist

clean.pip:
	pip uninstall locora --yes

clean.all:
	rm -rf *.egg-info
	rm -rf build
	rm -rf dist
	pip uninstall locora --yes
	find -name "*.pyc" -exec rm -f {} \;
