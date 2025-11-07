import copy

def print_matrix(matrix):
    is_permutation = all(all(elem == 0.0 or elem == 1.0 for elem in row) for row in matrix)

    for row in matrix:
        if is_permutation:
            print(" ".join(f"{int(elem):3d}" for elem in row))
        else:
            print(" ".join(f"{elem:9.4f}" for elem in row))
    print()


def lu_decomposition_pivoting(matrix_a, size):
    U = copy.deepcopy(matrix_a)
    L = [[0.0] * size for _ in range(size)]
    P = [[0.0] * size for _ in range(size)]

    for i in range(size):
        L[i][i] = 1.0
        P[i][i] = 1.0

    epsilon = 1e-9

    for j in range(size):
        pivot_row = j
        max_val = abs(U[j][j])

        for i in range(j + 1, size):
            if abs(U[i][j]) > max_val:
                max_val = abs(U[i][j])
                pivot_row = i

        if max_val < epsilon:
            raise ValueError("Matrix is singular.")

        if pivot_row != j:
            U[j], U[pivot_row] = U[pivot_row], U[j]
            P[j], P[pivot_row] = P[pivot_row], P[j]
            L[j][:j], L[pivot_row][:j] = L[pivot_row][:j], L[j][:j]

        for i in range(j + 1, size):
            factor = U[i][j] / U[j][j]
            L[i][j] = factor

            for k in range(j, size):
                U[i][k] -= (factor * U[j][k])

    return P, L, U


def main():
    size = 0
    while size <= 0:
        try:
            size = int(input("Enter the size of the square matrix (n x n): "))
            if size <= 0:
                print("Please enter a positive integer.")
        except ValueError:
            print("Invalid Value")

    a = [[0.0] * size for _ in range(size)]
    print(f"Enter your {size}x{size} matrix, one row at a time (numbers separated by spaces):")

    for i in range(size):
        while True:
            try:
                a[i] = list(map(float, input().strip().split()))
                if len(a[i]) != size:
                    print(f"Error: Please enter exactly {size} numbers.")
                    continue
                break
            except ValueError:
                print("Invalid input. Please enter numbers separated by spaces.")

    try:
        P, L, U = lu_decomposition_pivoting(a, size)

        print("\n--- P Matrix (Permutation) ---")
        print_matrix(P)

        print("--- L Matrix (Lower Triangular) ---")
        print_matrix(L)

        print("--- U Matrix (Upper Triangular) ---")
        print_matrix(U)

    except ValueError as e:
        print(f"\nError: {e}")


if __name__ == "__main__":
    main()