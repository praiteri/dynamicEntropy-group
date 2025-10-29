import logging

from ..orca_tools.extract_from_orca_output import extract_from_orca_output


class get_values(extract_from_orca_output):
    def read(self, **kwargs):
        self.line_indices = self.search_string(self.string)
        self.read_single_line()
        self.text = None
        if kwargs.get("text", None) is not None:
            self.text = kwargs["text"]

    def print_data(self, ll):
        separator = None
        if "..." in ll:
            ll = ll.replace("....", "...")
            separator = "..."
        elif ":" in ll:
            separator = ":"

        if separator is None:
            self.logger.info(ll)
        else:
            s = [x.strip() for x in ll.split(separator)]
            if self.text is None:
                self.logger.info(self.fmt.format(s[0], s[1]))
            else:
                self.logger.info(self.fmt.format(self.text, s[1]))

    def write(self, **kwargs):
        if self.result is None:
            return

        if len(self.line_indices) == 0:
            return

        if "first" in kwargs and kwargs["first"]:
            self.print_data(f"{self.result[0]}")

        elif "last" in kwargs and kwargs["last"]:
            self.print_data(f"{self.result[-1]}")

        else:
            for i in range(len(self.line_indices)):
                self.print_data(f"{self.result[i]}")
