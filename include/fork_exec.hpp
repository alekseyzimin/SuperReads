/**
 * Fork and exec, like execvp. Return the pid of the child process in
 * case of success or -1 in case of failure (of either fork, execvp or
 * some other syscall: close, pipe). Errno is set properly.
 */
pid_t fork_execvp(const char* file, char *const argv[]);

/**
 * Fork and exec, like execlp. See @fork_execvp.
 */
pid_t fork_execlp(const char* file, const char* arg, ...);

/**
 * Split command on white spaces, then fork and exec. See
 * @fork_execvp. Warning: this command modifies its input argument.
 */
pid_t fork_exec(char* cmd);

/**
 * Split command on white spaces, then fork and exec. See @fork_execvp.
 */
pid_t fork_exec(const char* cmd);
