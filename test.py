
def solution(n):

    dp = [[0 for _ in range(n+2)] for _ in range(n+2)]
    return recur(1, n, dp) - 1  # -1 the solution for one step with n brick.


def recur(height, target, dp):
    if dp[height][target] != 0:
        return dp[height][target]

    if target == 0:
        return 1
    if height > target:
        return 0

    dp[height][target] = recur(height+1, target, dp) + \
        recur(height+1, target-height, dp)
    return dp[height][target]


if __name__ == '__main__':
    print(solution(200))
