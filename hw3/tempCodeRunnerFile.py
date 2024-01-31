    print(k)
        v = matrix_vector_multiply(A, p)
        alpha_denominator = sum(a * b for a, b in zip(p, v))

        # Handle the case where the alpha denominator is close to zero
        if abs(alpha_denominator) < 1e-10:
            alpha = 0.0
        else:
            alpha = sum(a * b for a, b in zip(r, p)) / alpha_denominator

        x = [a + alpha * b for a, b in zip(x, p)]
        r_new = [a - alpha * b for a, b in zip(r, v)]
        z_new = solve_system(P, r_new)

        # Calculate beta denominator
        beta_denominator = sum(a * b for a, b in zip(r_new, z_new))

        # Handle the case where the beta denominator is close to zero
        if abs(beta_denominator) < 1e-10:
            beta = 0.0
        else:
            beta = beta_denominator / sum(a * b for a, b in zip(r, z))

        p = [a + beta * b for a, b in zip(z_new, p)]
        r = r_new
        z = z_new
        k += 1