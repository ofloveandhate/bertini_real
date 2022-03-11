from setuptools import find_packages, setup

EXCLUDE_FROM_PACKAGES = []

extras = {

    'optional': [
    ]
}

setup(name='bertini_real',
      version='1.6',  # TODO make this set programmatically
      description='Python library for bertini_real',
      url='http://bertinireal.com',
      author='Silviana Amethyst, with students Caden Joergens, Dan Hessler, Foong Min Wong',
      author_email='amethyst@uwec.edu',
      packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
      install_requires=['matplotlib',
                        'trimesh',
                        'pyopengl',
                        'glfw',
                        'triangle',
                        'dill',
                        'glumpy'],
      extras_require=extras,
      package_dir={'bertini_real': 'bertini_real'},
      zip_safe=False)

# on macos, i had to patch site-packages/OpenGL/platfomr/ctypesload.py
# as in https://githubmemory.com/repo/swistakm/pyimgui/issues/212
# that post indicated to change line 80, but I had to change line 35.
# it did work...  but it makes me uncomfortable