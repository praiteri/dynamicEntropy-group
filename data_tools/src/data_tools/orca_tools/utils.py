def is_numeric(fields):
    for x in fields:
        x = x.lower().strip()
        # Handle inf and nan explicitly
        if x not in ("inf", "-inf", "nan", "+inf", "infinity", "-infinity"):
            try:
                float(x)
            except ValueError:
                return False
        return True


def list_to_kwargs(input):
    """
    Converts a dictionary with a key that starts with "+" to a nested dictionary format.

    Input format: {'+animate': True, 'mode': 2, 'nframes': 3}
    Output format: {"animate": {'mode': 2, 'nframes': 3}}
    """
    result = {}
    nested_data = {}

    if type(input) is dict:
        for key, value in input.items():
            if key[0] == "+":
                result[key.replace("+", "")] = nested_data
            else:
                nested_data[key] = value
    elif type(input) is list:
        for ll in input:
            if ll[0] == "+":
                result[ll.replace("+", "")] = nested_data
            else:
                try:
                    k, v = ll.split("=")
                except ValueError:
                    k, v = ll, True
                nested_data[k] = v
    else:
        raise ValueError("input must be a list or a dictionary")
    return result


def clean_string(input_string):
    """
    Remove leading whitespace from each line in a string.

    Args:
        input_string: The input string to clean

    Returns:
        A cleaned string with leading characters removed from each line
    """
    cleaned_lines = []

    # Split the string into lines
    lines = input_string.strip().split("\n")

    for line in lines:
        # Remove leading whitespace first
        trimmed = line.lstrip()
        cleaned_lines.append(trimmed)

    # Join the cleaned lines back into a string
    return "\n".join(cleaned_lines) + "\n"


def centered_printout(text, width=40):
    parts = text.split("->")
    if len(parts) == 2:
        left, right = parts
        arrow = "->"
        formatted = f"{left:<8s}{arrow}{right:<8s}"
    else:
        # If there's no "->" in the string, center the whole text
        formatted = text.center(width)
    return formatted
