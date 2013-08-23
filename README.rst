BESL
====
:Author: `Brian Svoboda`_
:Email: svobodb@email.arizona.edu
:Source: https://github.com/autocorr/besl
:Docs: https://besl.readthedocs.org/en/latest
:Version: alpha

``besl`` is a general purpose library written in python for astronomical
research.


INSTALLING
----------
``besl`` is in alpha state so to install, you can clone it from the GitHub repository:

.. code-block::

    $ git clone git@github.com:autocorr/besl.git
    $ cd besl
    $ python setup.py install

or with ``pip``:

.. code-block::

    $ pip install git+http://github.com/autocorr/besl.git#egg=besl


REQUIREMENTS
------------
.. code-block::

    numpy      >= 1.7.1
    scipy      >= 0.10.1
    pandas     >= 0.11.0
    matplotlib >= 1.2.1
    astropy    >= 0.3.dev5057
    photutils  >= 0.0.dev128
    atpy       >= 0.9.6
    pyephem    >= 3.7.5.1
    pywcs      >= 1.11-4.8.2
    pyfits     >= 3.1.2


USING BESL
----------
Many paths in ``besl`` are hard-coded for specific directories and will need to be modified on the target machine. You can find info on how to use the sub-modules in the documentation at `ReadTheDocs`_.

To report a bug or request a feature, please use the issue tracker on `GitHub`_. Code contributions are welcome, just open a pull request.


LICENSE
-------
Copyright 2012, 2013 Brian Svoboda

BESL is free software: you can redistribute it and/or modify it under the terms
of the GNU General Public License (v3) as published by the Free Software
Foudnation, either version 3 of the License, or (at your option) any later
version.

BESL is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
BESL. If not, see http://www.gnu.org/licenses/.

.. _Brian Svoboda: http://autocorr.github.io
.. _ReadTheDocs: https://besl.readthedocs.org/en/latest
.. _GitHub: https://github.com/autocorr/besl
