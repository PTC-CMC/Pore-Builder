from setuptools import setup

setup(name='porebuilder',
        version='0.1',
        description='Graphene slit-pore builder for mBuild',
        url='http://github.com/rmatsum836/Pore-Builder',
        author='Ray Matsumoto',
        author_email='raymatsum@gmail.com',
        license='MIT',
        packages=['porebuilder'],
        zip_safe=False,
        entry_points={
            'mbuild.plugins':[
            "GraphenePore = porebuilder.porebuilder:GraphenePore",
            "GraphenePoreSolvent = porebuilder.porebuilder:GraphenePoreSolvent",
            ]
        }
    )
