import logging
from abc import ABC, abstractmethod
from ..logger import formatting


class extract_from_orca_output(ABC):
    def __init__(self, lines=None, string=None, **kwargs):
        # Set up the logger and
        # the base format for the logs
        self.logger = logging.getLogger("mylogger")
        self.mf = formatting()

        self.fmt = self.mf.fmt

        # All lines from the ORCA output file
        self.lines = lines
        self.line_indices = None
        if lines is None:
            raise Exception("No lines provided!")

        # Optional arguments
        self.kwargs = kwargs

        # String to search for
        self.string = string

        # Read the ORCA output and extract the data
        self.result = None
        self.read(**kwargs)

    def search_string(self, string, lines=None):
        """
        Search for a string in the lines
        """
        if lines is None:
            lines = self.lines
        # Extract line indices that match the string
        if len(string.split()) == 1:
            x = [i for i, line in enumerate(lines) if string in line.split()]
        else:
            x = [i for i, line in enumerate(lines) if string in line]
        if len(x) == 0:
            return None
        else:
            return x

    def parse_single_line(self, line_indices=None, lines=None):
        """
        Read a single line from the ORCA output file
        """
        if line_indices is None:
            if self.line_indices is None:
                return None
            line_indices = self.line_indices

        if lines is None:
            lines = self.lines

        # Initialize result
        result = None
        if len(line_indices) > 0:
            result = [lines[idx] for idx in line_indices]
        return result

    def read_single_line(self, line_indices=None, lines=None):
        self.result = self.parse_single_line(line_indices, lines)
        return

    def parse_blocks(self, string, escape_string, offset=0, line_indices=None, lines=None):
        """
        Read blocks of data from the ORCA output file
        """
        if lines is None:
            lines = self.lines

        if line_indices is None:
            line_indices = self.search_string(string, lines)

        if line_indices is None or len(line_indices) == 0:
            return None

        result = []
        for i in line_indices:
            idx = i + offset
            block = []

            def escape(string, elements):
                if not isinstance(elements, list):
                    elements = [elements]
                if any(element == "" for element in elements) and not string:
                    return True
                return any(element != "" and element in string for element in elements)

            while not escape(lines[idx], escape_string):
                block.append(lines[idx])
                idx += 1
            result.append(block)
        return result

    def read_blocks(self, string, escape_string, offset=0, line_indices=None, lines=None):
        self.result = self.parse_blocks(string, escape_string, offset, line_indices)
        return

    def result(self):
        """
        return the result of the search
        """
        return self.result

    @abstractmethod
    def read(self):
        pass

    @abstractmethod
    def write(self, **kwargs):
        pass
