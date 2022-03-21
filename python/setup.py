from setuptools import find_packages, setup

EXCLUDE_FROM_PACKAGES = []

extras = {
    'optional': [
    ]
}

about = {}
with open("bertini_real/__about__.py") as fp:
      exec(fp.read(), about)

setup(name='bertini_real',
      version=about['__version__'],  
      description=about['__summary__'],
      url=about['__uri__'],
      author=about['__author__'],
      author_email=about['__email__'],
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
      package_data={'bertini_real': ['surface/scad/*.scad']}, # for plugs and sockets on pieces of surfaces
      include_package_data=True, # for plugs and sockets on pieces of surfaces
      zip_safe=False)

# on macos, i had to patch site-packages/OpenGL/platfomr/ctypesload.py
# as in https://githubmemory.com/repo/swistakm/pyimgui/issues/212
# that post indicated to change line 80, but I had to change line 35.
# it did work...  but it makes me uncomfortable