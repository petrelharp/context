Troubleshooting
===============

1.  Non-finite likelihood values?

    Maybe the starting point isn't good,
    or maybe the scale on the parameters (e.g., `mutrates` or `selcoef`) is too large
    so that `optim` is trying to move too far away from your good starting point.
