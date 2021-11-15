from setuptools import find_packages, setup

setup(name='insilico_ab_lib',
        version='1.0',
        description='A script to generate in silico antibody sequence variants',
        url='--',
        author='MIT Lincoln Laboratory,
        license='BSD-2-Clause',
        packages=find_packages(),
        entry_points = {
            'console_scripts':['run_AbGen=scripts.run_library_gen:command_line'],
            },
        install_requires = [
            'argparse',
            'biopython',
            'numpy'
            ],
        zip_safe=True)

