// sat_solver.c
// Minimal DPLL-style SAT solver in C (expression tree based)

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Enum to represent types of expressions
typedef enum { VAR, CONST, AND, OR, NOT } ExprType;

// Expression node definition (expression tree)

typedef struct Expr {
    ExprType type; // Type of the expression node

    union {
        struct { struct Expr* left; struct Expr* right; } bin; // for AND, OR
        struct Expr* unary;   // for NOT
        int var_id;           // for variable nodes (e.g., a = 0, b = 1)
        int value;            // for constant nodes (0 or 1)
    } data;
} Expr;

// Forward declarations for functions
Expr* make_var(int var);
Expr* make_const(int value);
Expr* make_and(Expr* left, Expr* right);
Expr* make_or(Expr* left, Expr* right);
Expr* make_not(Expr* e);
Expr* simplify(Expr* e);
Expr* substitute(Expr* e, int var, int val);
int get_free_variable(Expr* e);
int is_satisfiable(Expr* e);
void free_expr(Expr* e);
Expr* parse_expr(const char* line);


// Constructors: allocate and return new nodes
Expr* make_var(int var) {
    Expr* e = malloc(sizeof(Expr));
    e->type = VAR;
    e->data.var_id = var;
    return e;
}

Expr* make_const(int value) {
    Expr* e = malloc(sizeof(Expr));
    e->type = CONST;
    e->data.value = value;
    return e;
}

Expr* make_and(Expr* left, Expr* right) {
    Expr* e = malloc(sizeof(Expr));
    e->type = AND;
    e->data.bin.left = left;
    e->data.bin.right = right;
    return e;
}

Expr* make_or(Expr* left, Expr* right) {
    Expr* e = malloc(sizeof(Expr));
    e->type = OR;
    e->data.bin.left = left;
    e->data.bin.right = right;
    return e;
}

Expr* make_not(Expr* inner) {
    Expr* e = malloc(sizeof(Expr));
    e->type = NOT;
    e->data.unary = inner;
    return e;
}


// Simplifies boolean expressions where possible
Expr* simplify(Expr* e) {
    if (e->type == CONST || e->type == VAR) return e;

    if (e->type == NOT) {
        Expr* inner = simplify(e->data.unary);
        if (inner->type == CONST) {
            int val = inner->data.value;
            free_expr(e);
            return make_const(!val);
        }
        e->data.unary = inner;
        return e;
    }

    Expr* l = simplify(e->data.bin.left);
    Expr* r = simplify(e->data.bin.right);

    if (e->type == AND) {

        // Apply basic simplification rules for AND
        if (l->type == CONST && l->data.value == 0) return make_const(0);
        if (r->type == CONST && r->data.value == 0) return make_const(0);
        if (l->type == CONST && l->data.value == 1) return r;
        if (r->type == CONST && r->data.value == 1) return l;
        e->data.bin.left = l;
        e->data.bin.right = r;
        return e;
    }

    if (e->type == OR) {
        
        // Apply simplification rules for OR
        if (l->type == CONST && l->data.value == 1) return make_const(1);
        if (r->type == CONST && r->data.value == 1) return make_const(1);
        if (l->type == CONST && l->data.value == 0) return r;
        if (r->type == CONST && r->data.value == 0) return l;
        e->data.bin.left = l;
        e->data.bin.right = r;
        return e;
    }

    return e;
}


// Replace variable `var` with constant `val` in expr

Expr* substitute(Expr* e, int var, int val) {
    if (e->type == VAR) {
        return (e->data.var_id == var) ? make_const(val) : make_var(e->data.var_id);
    }
    if (e->type == CONST) return make_const(e->data.value);
    if (e->type == NOT) return make_not(substitute(e->data.unary, var, val));
    return (e->type == AND) ?
        make_and(substitute(e->data.bin.left, var, val), substitute(e->data.bin.right, var, val)) :
        make_or(substitute(e->data.bin.left, var, val), substitute(e->data.bin.right, var, val));
}


// Returns the first variable found in the expr

int get_free_variable(Expr* e) {
    if (e->type == VAR) return e->data.var_id;
    if (e->type == CONST) return -1;
    if (e->type == NOT) return get_free_variable(e->data.unary);
    int l = get_free_variable(e->data.bin.left);
    int r = get_free_variable(e->data.bin.right);
    return (l != -1) ? l : r;
}


// DPLL-style recursive SAT checking

int is_satisfiable(Expr* e) {
    e = simplify(e); // simplify before checking
    if (e->type == CONST) return e->data.value; // base case: true/false

    int var = get_free_variable(e);
    if (var == -1) return 0; // no variables left but not a constant => contradiction

    // Try assigning var = 1 and var = 0
    Expr* try_true = substitute(e, var, 1);
    Expr* try_false = substitute(e, var, 0);

    int result = is_satisfiable(try_true) || is_satisfiable(try_false);

    free_expr(try_true);
    free_expr(try_false);
    return result;
}


// Recursively free expression tree

void free_expr(Expr* e) {
    if (!e) return;
    if (e->type == AND || e->type == OR) {
        free_expr(e->data.bin.left);
        free_expr(e->data.bin.right);
    } else if (e->type == NOT) {
        free_expr(e->data.unary);
    }
    free(e);
}


// Parse postfix expression string (e.g. "a b &")
// Supports: a-z, 0/1, &, |, !

Expr* parse_expr(const char* line) {
    const char* p = line;
    Expr* stack[256];
    int sp = 0; // stack pointer

    while (*p) {
        if (*p == ' ') { p++; continue; }

        // Push variable (a to z => var_id 0–25)
        else if (*p >= 'a' && *p <= 'z') {
            stack[sp++] = make_var(*p - 'a');
        }

        // Push constant 0 or 1
        else if (*p == '0' || *p == '1') {
            stack[sp++] = make_const(*p - '0');
        }

        // Unary NOT: !a
        else if (*p == '!') {
            p++;
            if (*p >= 'a' && *p <= 'z') {
                stack[sp++] = make_not(make_var(*p - 'a'));
            }
        }

        // Binary AND
        else if (*p == '&' && sp >= 2) {
            Expr* r = stack[--sp];
            Expr* l = stack[--sp];
            stack[sp++] = make_and(l, r);
        }

        // Binary OR
        else if (*p == '|' && sp >= 2) {
            Expr* r = stack[--sp];
            Expr* l = stack[--sp];
            stack[sp++] = make_or(l, r);
        }

        p++;
    }

    // At the end, there should be exactly one expression on the stack
    return (sp == 1) ? stack[0] : NULL;
}

// MAIN FUNCTION — reads expressions from file
// and prints SAT/UNSAT for each
int main() {
    FILE* f = fopen("expression.txt", "r");
    if (!f) {
        perror("expression.txt");
        return 1;
    }

    char line[256];
    int count = 1;
    while (fgets(line, sizeof(line), f)) {
        line[strcspn(line, "\n")] = 0; // Remove newline

        // Strip comments (start from "//")
        char* comment = strstr(line, "//");
        if (comment) *comment = '\0';

        // Trim trailing spaces
        for (int i = strlen(line) - 1; i >= 0 && line[i] == ' '; i--) {
            line[i] = '\0';
        }

        if (strlen(line) == 0) continue;

        Expr* e = parse_expr(line);
        if (!e) {
            printf("[%d] Parse error: %s\n", count, line);
        } else {
            printf("[%d] %s => %s\n", count, line, is_satisfiable(e) ? "SAT" : "UNSAT");
            free_expr(e);
        }
        count++;
    }

    fclose(f);
    return 0;
}