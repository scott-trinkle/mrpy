from setuptools import setup


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='mrpy',
      version='0.1',
      description='Python package wrapping MRTrix functions',
      long_description=readme(),
      author='Scott Trinkle',
      author_email='tscott.trinkle@gmail.com',
      license='MIT',
      packages=['mrpy'],
      package_dir={'mrpy': 'mrpy'},
      package_data={'mrpy': ['data/*']},
      zip_safe=False)
