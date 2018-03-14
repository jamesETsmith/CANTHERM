# CANTHERM
Python based software to calculate thermodynamic properties of molecules and rate coefficients of reactions. At the current time Gaussian log files are require to for the hindered rotor scans and frequency calculations and the CanTherm can extract the energy from either Gaussian or MOLPRO files.

# Rules for Contributing
If you contribute to this repository, we are trying to follow the Google Python style guide as closely as possible so please make sure to use the naming conventions described below. If you'd like to know more about the style guide please see this [link](https://google.github.io/styleguide/pyguide.html).

| Type                       | Public             | Internal                                                          |
|:---------------------------|:------------------:|:-----------------------------------------------------------------:|
| Packages                   | lower_with_under   |                                                                   |
| Modules                    | lower_with_under   | _lower_with_under                                                 |
| Classes                    | CapWords           | _CapWords                                                         |
| Exceptions                 | CapWords           |                                                                   |
| Functions	                 | lower_with_under() | _lower_with_under()                                               |
| Global/Class Constants     | CAPS_WITH_UNDER    |	_CAPS_WITH_UNDER                                                  |
| Global/Class Variables     | lower_with_under   |	_lower_with_under                                                 |
| Instance Variables         | lower_with_under   |	_lower_with_under (protected) or __lower_with_under (private)     |
| Method Names               | lower_with_under() | _lower_with_under() (protected) or __lower_with_under() (private) |
| Function/Method Parameters | lower_with_under	  |                                                                   |
| Local Variables            | lower_with_under   |                                                                   |
