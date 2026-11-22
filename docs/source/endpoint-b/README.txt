Endpoint-B negative direction in Theorem 2.1 of the revision.

In a Python 3.11+ environment, install the version in requirements.txt.
Run python verify_endpoint_b.py from this directory. The checker is
independent of the historical project scripts and reads no external data.

Expected certificate result at 768 bits
  -1.119e-28 < quadratic form < -1.118e-28
  support strictly beyond the active trace interval [0.02, 0.545]

The coefficients in verify_endpoint_b.py are exact decimals. The checker
evaluates the defining full endpoint-B kernel via the incomplete gamma
moments in equations (2.5) and (2.6), with a rigorous remainder bound for
logarithmic moments. The manuscript proves the passage from node values
to the smooth closed-trace witness and the consequence for minimizers.
