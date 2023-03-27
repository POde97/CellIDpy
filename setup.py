from setuptools import setup

setup(
    name='Cell-IDpy',
    version='1.0',    
    description='mca based method for cell type labeling',
    url=https://github.com/POde97/Cell-IDpy,
    author='Paolo Odello',
    author_email='paoloodeo.o@gmail.com',
    license='MIT license',
    packages=['CellID'],
    install_requires=['mapply==0.1.21',
                      'scanpy==1.9.3',                     
                      ],
)
