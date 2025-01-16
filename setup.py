#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name="DiPPER2",  # The name of your package (as it will appear on PyPI)
    version="0.1.0",           # Initial version of your package
    author="Theresa Wacker",        
    author_email="thwa223@proton.me", 
    description="Diagnostic Primer Picking and Evaluation pipeline for Reliability and Reproducibility",  
    long_description=open('README.md').read(),  
    long_description_content_type="text/markdown",  
    url="https://github.com/ThWacker/DiPPER2",  
    packages=find_packages(),  
    classifiers=[  
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0-only",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',  
    install_requires=[  
            'bio>=1.6.2', 
            'nbconvert>=7.11.0,<8.0.0'
    ],
    tests_require=['unittest'],  
    test_suite='tests',  
    
)
