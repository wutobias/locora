from setuptools import setup, Extension, find_packages
import numpy as np

__version__ = "0.2"

# define the extension module
extensions = []
extensions.append(Extension('locora.grid_solvent._crd_systems_ext',
                            sources=['./src/_crd_systems_ext.c',
                                     './src/Vec.c'],
                            include_dirs=[np.get_include()],
                            language='c'))

extensions.append(Extension('locora.process_data._utils_ext',
                            sources=['./src/_process_data_ext.c',
                                     './src/Vec.c'],
                            include_dirs=[np.get_include()],
                            language='c'))

extensions.append(Extension('locora.utils._read_write_ext',
                            sources=['./src/_read_write_ext.c'],
                            include_dirs=[np.get_include()],
                            language='c'))

setup(name='locora',
      author='Tobias Wulsdorf',
      author_email='tobias.wulsdorf@gmail.com',
      description='LoCorA: A tool for studying local correlations between solute entities and solvent molecules',
      version=__version__,
      license='MIT',
      platforms=['Linux'],
      packages=find_packages(),
      ext_modules=extensions,
      zip_safe=False,
      entry_points={
          'console_scripts':
              ['run_locora  = locora.scripts.run_locora:entry_point']}, )
