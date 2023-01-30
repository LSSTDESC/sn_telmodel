from setuptools import setup

# get the version here
pkg_vars = {}

with open("version.py") as fp:
    exec(fp.read(), pkg_vars)

setup(
    name='sn_telmodel',
    version=pkg_vars['__version__'],
    description='Set of tools for telescope model',
    url='http://github.com/lsstdesc/sn_telmodel',
    author='Philippe Gris',
    author_email='philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages= ['sn_telmodel'],
    # All files from folder sn_script_input
    # package_data={'sn_script_input': ['*.txt']},
    python_requires='>=3.5',
    zip_safe=False
)
