from setuptools import setup

setup(
    name='CellIDpy',
    version='1.0',    
    description='mca based method for cell type labeling',
    url='https://github.com/POde97/CellIDpy',
    author='Paolo Odello',
    author_email='paoloodeo.o@gmail.com',
    package_data={'CellIDpy': ['HgProteinCoding', 'MegProteinCoding']}
    license='MIT license',
    packages=['CellIDpy'],
    install_requires=['mapply==0.1.21',
                      'scanpy==1.9.3',                     
                      ],
)
