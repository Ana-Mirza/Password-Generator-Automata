#include <bits/stdc++.h>

using namespace std;

#define KMAX 205
#define MOD  1009

int n;
int adj[KMAX][KMAX];

int odd_power, even_power;
bool multiplied_tmp = false;
bool is_even = false;
bool is_odd = false;
bool found_solution = false;

// intermediate matrices
vector<vector<vector<int>>> even_matrix = vector<vector<vector<int>>>();
vector<vector<vector<int>>> odd_matrix = vector<vector<vector<int>>>();
map<int, int> odd_matrix_index;

vector<map<int, string>> delta_function = vector<map<int, string>>(5001, map<int, string>());
map<int, string> alphabet = {{1, "a"}, {2, "b"}, {3, "c"}, {4, "d"}, {5, "e"},
                        {6, "f"}, {7, "g"}, {8, "h"}, {9, "i"}, {10, "j"},
                        {11, "k"}, {12, "l"}, {13, "m"}, {14, "n"}, {15, "o"},
                        {16, "p"}, {17, "q"}, {18, "r"}, {19, "s"}, {20, "t"},
                        {21, "u"}, {22, "v"}, {23, "w"}, {24, "x"}, {25, "y"},
                        {26, "z"}};

// C = A * B
vector<vector<int>> intermediates = vector<vector<int>>(205, vector<int>(205, 0));
void multiply_matrix(int A[KMAX][KMAX], int B[KMAX][KMAX], int C[KMAX][KMAX]) {
    int tmp[KMAX][KMAX];
    intermediates = vector<vector<int>>(205, vector<int>(205, 0));

    // tmp = A * B
    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= n; ++j) {
            unsigned long long sum = 0;
 
            for (int k = 1; k <= n; ++k) {
                sum += A[i][k] * B[k][j];
                if (sum > 100)
                    sum = 1;
                if (A[i][k] * B[k][j] != 0) {
                    intermediates[i][j] = k;
                }
            }

            if (sum != 0)
                tmp[i][j] = 1;
            else
                tmp[i][j] = 0;
        }
    }

    // copy last auxiliary matrix if it is the first odd power
    if (even_power != 0 && !multiplied_tmp && is_odd && found_solution) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                int k = even_matrix.size() - 1;
                intermediates[i][j] = even_matrix[k][i][j];
            }
        }
    }

    // copy adj matrix if there is no even power previous to it
    if (!multiplied_tmp && is_odd && even_power == 0 && found_solution) {
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                intermediates[i][j] = adj[i][j];
            }
        }
    }
 
    //  C = tmp
    memcpy(C, tmp, sizeof(tmp));
}
 
// R = C^p
void power_matrix(int C[KMAX][KMAX], int p, int R[KMAX][KMAX]) {
    // tmp = I (identity matrix)
    int tmp[KMAX][KMAX];
    for (int i = 0; i < KMAX; ++i) {
        for (int j = 0; j < KMAX; ++j) {
            tmp[i][j] = (i == j) ? 1 : 0;
        }
    }

    even_power = 0;
    odd_power = 0;
    multiplied_tmp = false;
 
    while (p != 1) {
        if  (p % 2 == 0) {
            // determine even power
            is_even = true;
            if (even_power == 0)
                even_power = 2;
            else 
                even_power *= 2;

            multiply_matrix(C, C, C);     // C = C*C
            p /= 2;                       // we have left C^(p/2) to compute
            is_even = false;

            // add matrix in vector
            if (found_solution)
                even_matrix.push_back(intermediates);
        } else {
            // determine odd power
            is_odd = true;
            if (even_power == 0 && odd_power == 0)
                odd_power = 1;
            else
                odd_power += even_power;

            // reduce to previous case:
            multiply_matrix(tmp, C, tmp); // tmp = tmp*C
            --p;                          // we have left C^(p-1) to compute
            multiplied_tmp = true;
            is_odd = false;

            // add matrix in vector
            if (found_solution) {
                odd_matrix_index[odd_power] = odd_matrix.size();
                odd_matrix.push_back(intermediates);
            }
        }
    }
 
    if (odd_power == 0)
        is_odd = true;
    // solution is in C and in tmp
    multiply_matrix(tmp, C, R);           // result = tmp * C
    is_odd = false;

    // add last matrix of intermediates
    if (found_solution) {
        odd_matrix_index[odd_power + even_power] = odd_matrix.size();
        odd_matrix.push_back(intermediates);
    }
}

string square_power_decomposition(int left_node, int right_node, int length) {
    // leaf: return tranzition symbol from tranzition function
    if (length == 1)
        return delta_function[left_node][right_node];

    // find left part and right part of the word with half the length and combine them
    int new_len = length / 2;
    int index = (int)log2(length) - 1;
    int node = even_matrix[index][left_node][right_node];
    return square_power_decomposition(left_node, node, new_len) + square_power_decomposition(node, right_node, new_len);
}

string power_decomposition(int left_node, int right_node, int length) {
    if (length == 0)
        return "";

    // leaf: return tranzition symbol from tranzition function
    if (length == 1)
        return delta_function[left_node][right_node];

    // find greatest power of 2 less than current lenght
    int right_len = exp2((int)log2(length));
    int left_len = length - right_len;
    // check if length was a power of 2
    if (left_len == 0) {
        right_len /= 2;
        left_len = right_len;
    }

    int intermediate_node = odd_matrix[odd_matrix_index[length]][left_node][right_node];

    string left_word;
    // find left part of the word
    if ((left_len & (left_len - 1)) == 0)
        left_word = square_power_decomposition(left_node, intermediate_node, left_len);
    else
        left_word = power_decomposition(left_node, intermediate_node, left_len);

    // find right part of the word
    string right_word = square_power_decomposition(intermediate_node, right_node, right_len);

    // combine solutions
    return left_word + right_word;
}

void task2() {
    freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);

    int t, sigma, m, k;
    int init_node;
    unordered_set<int> fin;
    delta_function = vector<map<int, string>>(5001, map<int, string>());

    // read parameters
    cin >> t >> n >> sigma >> m >> k;
    // read initial state
    cin >> init_node;

    // read final states
    for (int i = 0; i < m; i++) {
        int val;
        cin >> val;
        fin.insert(val);
    }

    // read transitions
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= sigma; j++) {
            int val;
            cin >> val;
            adj[i][val] = 1;
            delta_function[i][val] = alphabet[j];
        }
    }

    // find minimum length greater than or equal to k
    int k_tmp = k - 1;
    int final_node, len;

    while (k_tmp < k + n && !found_solution) {
        k_tmp++;
        int R[KMAX][KMAX];
        int tmp[KMAX][KMAX];
        memcpy(tmp, adj, sizeof(adj));

        // find paths of length k_tmp
        power_matrix(tmp, k_tmp, R);

        // check if found solution
        for (int i = 1; i <= n; i++) {
            if (R[init_node][i] != 0 && fin.find(i) != fin.end() && k_tmp >= k) {
                final_node = i;
                len = k_tmp;
                cout << k_tmp << "\n";
                // solution = true;
                found_solution = true;
                break;
            }
        }
    }

    // solution not found
    if (!found_solution) {
        cout << -1 << "\n";
        return;
    }

    // find intermediate nodes
    int R[KMAX][KMAX];
    int tmp[KMAX][KMAX];
    memcpy(tmp, adj, sizeof(adj));
    power_matrix(tmp, len, R);

    // reconstruct word
    string word = power_decomposition(init_node, final_node, len);
    cout << word << "\n";
    return;
}

void task1() {
    freopen("input.txt", "r", stdin);
	freopen("output.txt", "w", stdout);
    ios_base::sync_with_stdio(false);

    int t, n, sigma, m, k;
    int init;
    unordered_set<int> fin;
    vector<vector<int>> dp = vector<vector<int>>(25002, vector<int>(5001, 0));
    vector<vector<int>> neighbor = vector<vector<int>>(5001, vector<int>());

    // read parameters
    cin >> t >> n >> sigma >> m >> k;
    // read initial state
    cin >> init;

    // read final states
    for (int i = 0; i < m; i++) {
        int val;
        cin >> val;
        fin.insert(val);
    }

    // read transitions
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= sigma; j++) {
            int val;
            cin >> val;
            neighbor[i].push_back(val);
            delta_function[i][val] = alphabet[j];
        }
    }

    // find word of minimum length equal to at least k
    int k_tmp = 1, iterations = 0;
    bool solution = false;
    int end_state;

    // dp base case
    for (auto neigh : neighbor[init]) {
        dp[1][neigh] = init;
    }

    // complete dp matrix
    while (k_tmp < k || !solution) {
        k_tmp++;

        // count interations after minimum length
        if (k_tmp >= k) {
            iterations++;
        }

        for (int i = 1; i <= n; i++) {
            // extend solution
            if (dp[k_tmp - 1][i] != 0) {
                for (auto neigh : neighbor[i]) {
                    if (dp[k_tmp][neigh] == 0)
                        dp[k_tmp][neigh] = i;

                    if (k_tmp >= k && fin.find(neigh) != fin.end()) {
                        solution = true;
                        end_state = neigh;
                        break;
                    }
                }
            }

            // check for cycle
            // if there is no final solution in n steps, we have a cycle
            // maximum length between two nodes is n - 1
            if (iterations > n) {
                cout << -1 << "\n";
                return;
            }
        }
    }

    // reconstruct word by going backwards from parent to parent
    cout << k_tmp << "\n";
    string word = "";
    while (k_tmp > 0) {
        int parent = dp[k_tmp][end_state];
        word += delta_function[parent][end_state];
        end_state = parent;
        k_tmp--;
    }
    reverse(word.begin(), word.end());
    cout << word << "\n";
}

int main() {
    freopen("input.txt", "r", stdin);

    int t;
    cin >> t;

    if (t == 1) {
        task1();
    } else if (t == 2) {
        task2();
    }

    return 0;
}