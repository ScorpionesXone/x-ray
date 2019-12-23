from setuptools import setup, find_packages

with open('READ.txt') as fh:
    long_description = fh.read()

setup(
    name='Xray4U',
    version='0.0.1',
    author='Arthur L. Avakyan',
    author_email='artur_ava97@mail.ru',
    description='Light Curve Approximation with thermal wind',
    long_description=long_description,
    long_description_content_type='text/markdown',
    #long_description=open(join(dirname(__file__), 'READ.txt')).read(),
    url='http://phys.msu.ru',
    license='SAI',
    packages=find_packages(),
    #package_data = ['data/asu.tsv'],
    #scripts=['module/rytp.py'],
    test_suite='test',
    entry_points={
        'console_scripts': ['rytp = module.rytp:main'],
        #'gui_scripts': ['plotsnake = ser.snake:plot'],
    },
    install_requires=['numpy>=1.13', 'scipy>=1.0', 'matplotlib>=2.0,<3.0.3'],
    keywords='accretion, accretion disks – binaries, thermal winds, black holes, outflows – X-Ray : binaries',
)
