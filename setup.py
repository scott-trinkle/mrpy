from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='mrpy',
      version='0.1',
      description='Python package wrapping MRTrix functions',
      long_description=readme(),
      author='Chineze Egwudo and Scott Trinkle',
      author_email='tscott.trinkle@gmail.com',
      license='MIT',
      packages=['mrpy'],
      zip_safe=False)
