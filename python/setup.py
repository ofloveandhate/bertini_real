from setuptools import find_packages, setup

EXCLUDE_FROM_PACKAGES = []

extras = {

    'optional': [
    ]
}

setup(name='bertini_real',
      version='1.7',  # TODO make this set programmatically
      description='Python library for bertini_real',
      url='https://bertinireal.com',
      author='silviana amethyst, with students Caden Joergens, Dan Hessler, Foong Min Wong',
      author_email='amethyst@uwec.edu',
      packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
      install_requires=['matplotlib',
                        'trimesh',
                        'pyopengl',
                        'glfw',
                        'triangle',
                        'dill',
                        'glumpy',
                        'scipy'],
      extras_require=extras,
      package_dir={'bertini_real': 'bertini_real'},
      package_data={'bertini_real': ['surface/scad/*.scad']},
      include_package_data=True,
      zip_safe=False)

# on macos, i had to patch site-packages/OpenGL/platfomr/ctypesload.py
# as in https://githubmemory.com/repo/swistakm/pyimgui/issues/212
# that post indicated to change line 80, but I had to change line 35.
# it did work...  but it makes me uncomfortable