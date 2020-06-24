#!/usr/bin/env python
"""
********************************************************************************
*                                                                              *
*    CANTHERM                                                                  *
*                                                                              *
*    Copyright (C) 2020                                                        *
*    Sandeep Sharma, James E. T. Smith, and Virginia Johnson                   *
*                                                                              *
*    This program is free software: you can redistribute it and/or modify      *
*    it under the terms of the GNU General Public License as published by      *
*    the Free Software Foundation, either version 3 of the License, or         *
*    (at your option) any later version.                                       *
*                                                                              *
*    This program is distributed in the hope that it will be useful,           *
*    but WITHOUT ANY WARRANTY; without even the implied warranty of            *
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
*    GNU General Public License for more details.                              *
*                                                                              *
*    You should have received a copy of the GNU General Public License         *
*    along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
*                                                                              *
********************************************************************************
"""
import os


def get_sample_file_path(filename: str) -> str:
    """Gets the full file path for any sample file in the data directory.

    Parameters
    ----------
    filename : str
        The relative file path of the desire sample file. E.g. "bz.log".

    Returns
    -------
    str
        The absolute file path for the sample file requested.

    Raises
    ------
    ValueError
        If the filename doesn't exist in the data directory.
    """

    cantherm_path = os.path.abspath(os.path.dirname(__file__))
    full_path = os.path.join(cantherm_path, "..", "data", filename)

    if not os.path.exists(full_path):
        raise ValueError(
            f"The filename you've chosen ({filename}) doesn't exist in the Cantherm sample files"
        )

    return full_path

