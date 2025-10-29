def format_output_dict(
    statistics_dict, key_width=25, decimal_places=6, sci_threshold=1e6
):
    """
    Format a dictionary with right-aligned keys and decimal-aligned numeric values.
    Handles both regular notation and scientific notation.

    Args:
        statistics_dict (dict): Dictionary to format
        key_width (int): Width for the key column (default: 25)
        decimal_places (int): Number of decimal places to display (default: 6)
        sci_threshold (float): Threshold for switching to scientific notation (default: 1e6)

    Returns:
        str: Formatted string representation of the dictionary
    """
    formatted_lines = []

    # Convert all values to strings with appropriate formatting
    formatted_values = []
    for key, value in statistics_dict.items():
        if isinstance(value, (int, float)):
            # Decide whether to use scientific notation
            abs_value = abs(float(value))
            if abs_value >= sci_threshold or (
                abs_value != 0 and abs_value < 1 / sci_threshold
            ):
                # Use scientific notation
                formatted_value = f"{float(value):.{decimal_places}e}"
            else:
                # Use regular notation
                formatted_value = f"{float(value):.{decimal_places}f}"
        else:
            # Handle non-numeric values
            formatted_value = str(value)
        formatted_values.append((key, formatted_value))

    # Separate regular numbers from scientific notation
    regular_values = []
    sci_values = []

    for key, formatted_value in formatted_values:
        if "e" in formatted_value.lower():
            sci_values.append((key, formatted_value))
        else:
            regular_values.append((key, formatted_value))

    # Find alignment for regular numbers (decimal point alignment)
    max_before_decimal = 0
    if regular_values:
        for _, formatted_value in regular_values:
            if "." in formatted_value:
                before_decimal = len(formatted_value.split(".")[0])
                max_before_decimal = max(max_before_decimal, before_decimal)

    # Find alignment for scientific notation (align the 'e')
    max_before_e = 0
    if sci_values:
        for _, formatted_value in sci_values:
            if "e" in formatted_value.lower():
                before_e = len(formatted_value.split("e")[0])
                max_before_e = max(max_before_e, before_e)

    # Determine the overall value column width
    regular_width = 0
    sci_width = 0

    if regular_values:
        # Width needed for regular numbers
        sample_regular = regular_values[0][1]
        if "." in sample_regular:
            after_decimal = len(sample_regular.split(".")[1])
            regular_width = (
                max_before_decimal + 1 + after_decimal
            )  # +1 for decimal point
        else:
            regular_width = max_before_decimal

    if sci_values:
        # Width needed for scientific notation
        sample_sci = sci_values[0][1]
        after_e = len(sample_sci.split("e")[1])
        sci_width = max_before_e + 1 + after_e  # +1 for 'e'

    # Use the maximum width for consistent alignment
    value_width = max(regular_width, sci_width)

    # Format each line
    all_formatted = []

    # Process regular numbers
    for key, formatted_value in regular_values:
        aligned_key = key.rjust(key_width)

        if "." in formatted_value:
            before_decimal, after_decimal = formatted_value.split(".")
            padded_before = before_decimal.rjust(max_before_decimal)
            aligned_value = f"{padded_before}.{after_decimal}"
        else:
            aligned_value = formatted_value.rjust(max_before_decimal)

        # Right-align within the overall value width
        final_value = aligned_value.rjust(value_width)
        all_formatted.append((key, f"{aligned_key} : {final_value}"))

    # Process scientific notation numbers
    for key, formatted_value in sci_values:
        aligned_key = key.rjust(key_width)

        if "e" in formatted_value.lower():
            before_e, after_e = formatted_value.split("e")
            padded_before = before_e.rjust(max_before_e)
            aligned_value = f"{padded_before}e{after_e}"
        else:
            aligned_value = formatted_value

        # Right-align within the overall value width
        final_value = aligned_value.rjust(value_width)
        all_formatted.append((key, f"{aligned_key}: {final_value}"))

    # Sort back to original order
    result_dict = {item[0]: item[1] for item in all_formatted}
    formatted_lines = [result_dict[key] for key in statistics_dict.keys()]

    return formatted_lines
    # return "\n".join(formatted_lines)


# Alternative function with automatic threshold detection
def format_output_dict_auto(statistics_dict, key_width=25, decimal_places=6):
    """
    Automatically determine the best formatting based on the range of values.
    """
    # Analyze the values to determine appropriate formatting
    numeric_values = [
        v for v in statistics_dict.values() if isinstance(v, (int, float))
    ]

    if not numeric_values:
        return format_output_dict(statistics_dict, key_width, decimal_places, 1e6)

    max_abs = max(abs(v) for v in numeric_values if v != 0)
    min_abs = (
        min(abs(v) for v in numeric_values if v != 0)
        if any(v != 0 for v in numeric_values)
        else 1
    )

    # Use scientific notation if we have very large or very small numbers
    if max_abs >= 1e6 or min_abs <= 1e-4:
        threshold = 1e4  # Lower threshold to group more numbers in scientific notation
    else:
        threshold = 1e6  # Higher threshold to keep most numbers in regular notation

    return format_output_dict(statistics_dict, key_width, decimal_places, threshold)
