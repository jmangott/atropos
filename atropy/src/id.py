"""
Contains the `Id` class for proper naming of the partitions of the reaction network.
"""

import re
from typing import Union


class Id(str):
    def __init__(self, _id: str = ""):
        if re.match("^[01]*$", _id):
            self.id = _id
        else:
            raise ValueError("Not a valid input string")

    def __str__(self):
        if self.id == "":
            return "root"
        else:
            return self.id

    def __add__(self: "Id", other: Union["Id", int]):
        if isinstance(other, Id) and other.id in ("0", "1"):
            return Id(self.id + other.id)
        elif isinstance(other, int) and other in (0, 1):
            return Id(self.id + str(other))
        raise ValueError("Addition not defined")

    def __sub__(self: "Id", other: Union["Id", int]):
        if (isinstance(other, Id) and other.id == "1") or (
            isinstance(other, int) and other == 1
        ):
            if self.id == "":
                raise ValueError("Subtraction from 'root' not defined")
            else:
                return Id(self.id[:-1])
        else:
            raise ValueError("Subtraction not defined")

    def __int__(self: "Id"):
        return int("0b" + self.id, 2)


if __name__ == "__main__":
    a = Id("")
    print(str(a))
