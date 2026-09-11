def validate_radius(radius):
    if radius <= 5:
        raise ValueError("Radius must be > 5 A")


def validate_shell(shell, radius):
    if shell <= 0:
        raise ValueError("Shell thickness must be > 0")
    if shell >= radius:
        raise ValueError("Shell thickness must be smaller than radius")


def validate_disorder(disorder, shell):
    if disorder <= 0:
        raise ValueError("Disorder must be > 0")
    if disorder > shell:
        raise ValueError("Disorder should be smaller than shell thickness")