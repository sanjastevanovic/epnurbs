[metadata]
# This includes the license file in the wheel.
license_file = LICENSE.txt

[bdist_wheel]
# This flag says to generate wheels that support both Python 2 and Python
# 3. If your code will not run unchanged on both Python 2 and 3, you will
# need to generate separate wheels for each Python version that you
# support. Removing this line (or setting universal to 0) will prevent
# bdist_wheel from trying to make a universal wheel. For more see:
# https://packaging.python.org/tutorials/distributing-packages/#wheels
universal=0
