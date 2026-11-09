"""
Cheap consistency check between an entry's LaTeX equations and its declared symbols.

It catches renamed, misspelled, undeclared or unused symbols. It does NOT
check the mathematics: a wrong sign or operator passes. Each model's
test_validation.py is what checks that rhs is correct.
"""

import re

# LaTeX commands allowed in an equation besides the declared symbols.
ALLOWED_COMMANDS = {
    "frac", "left", "right", "cdot", "times", "exp", "ln", "log", "sqrt", "quad",
}
# Characters allowed besides symbols: digits, operators, grouping, whitespace.
ALLOWED_CHARS = set("0123456789.+-*/=()[]{}^_,| \t\n")
# Symbols every equation may use without declaring them.
BUILTIN_SYMBOLS = ("t",)


def _tokenizer(symbols):
    ordered = sorted(set(symbols) | set(BUILTIN_SYMBOLS), key=len, reverse=True)
    alternatives = [re.escape(symbol) for symbol in ordered]
    alternatives += [r"\\[A-Za-z]+", r"\\.", r"[A-Za-z]", r"."]
    return re.compile("|".join(f"(?:{a})" for a in alternatives), re.DOTALL)


def symbols_used(equation, symbols):
    """
    Return (used, problems): the declared symbols an equation uses, and every
    token that is neither a declared symbol nor allowed LaTeX.
    """
    declared = set(symbols) | set(BUILTIN_SYMBOLS)
    used, problems = set(), []
    for token in _tokenizer(symbols).findall(equation):
        if token in declared:
            used.add(token)
        elif token.startswith("\\") and len(token) > 1 and token[1:].isalpha():
            if token[1:] not in ALLOWED_COMMANDS:
                problems.append(f"unknown LaTeX command or undeclared symbol {token!r}")
        elif token.isalpha():
            problems.append(f"undeclared symbol {token!r}")
        elif token not in ALLOWED_CHARS:
            problems.append(f"unexpected character {token!r}")
    return used & set(symbols), problems
