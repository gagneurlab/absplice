from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()


requirements = [
    'setuptools',
]
extras_requirements = {
    "predict": [
        'kipoiseq>=0.3.0',
        'mmsplice>=2.1.0'
    ]
}

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest', 'pytest-benchmark']

setup(
    author="M. Hasan Celik, Nils Wagner",
    author_email='muhammedhasancelik@gmail.com, wagnern@in.tum.de',
    classifiers=[
        'Development Status :: 1 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="Aberrant splicing prediction across human tissues",
    install_requires=requirements,
    license="MIT license",
    entry_points={
        'console_scripts': [
            'spliceai_rocksdb_download=absplice.spliceai_rocksdb_download:spliceai_rocksdb_download'
        ]
    },
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/markdown',
    include_package_data=True,
    keywords='absplice',
    name='absplice',
    packages=find_packages(include=['absplice']),
    setup_requires=setup_requirements,
    extras_require=extras_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/gagneurlab/absplice',
    version='0.0.1',
    zip_safe=False
)
