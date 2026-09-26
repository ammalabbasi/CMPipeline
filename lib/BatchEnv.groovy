import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption

/**
 * The batch-correction conda lockfile, ready for Nextflow to build the environment.
 *
 * ConQuR (ivartb/ConQuR_par) and cqrReg 1.2.1 are on no channel the pipeline uses, so they are
 * built from conda_envs/recipes/ and committed under conda_envs/local_channel/. An explicit
 * lockfile needs absolute URLs; the committed lockfile names these two packages under
 * @LOCAL_CHANNEL@, replaced here by file://<projectDir>/conda_envs/local_channel. The result is
 * written to <workDir>/conda_locks/ (the checkout may be read-only), rewritten only when its
 * content changes, so -resume and the conda cache see a stable file.
 */
class BatchEnv {
    static final String PLACEHOLDER = '@LOCAL_CHANNEL@'
    static final String LOCKFILE = 'conda_envs/lock/batch_correction_env.linux-64.txt'

    static String lockfile(Object projectDir, Object workDir) {
        Path root = (projectDir as Path).toAbsolutePath().normalize()
        String template = root.resolve(LOCKFILE).text
        if (!template.contains(PLACEHOLDER))
            throw new IllegalStateException("${LOCKFILE} does not reference the local channel")
        String channel = root.resolve('conda_envs/local_channel').toUri().toString().replaceAll('/+$', '')
        String text = template.replace(PLACEHOLDER, channel)
        Path dir = (workDir as Path).toAbsolutePath().normalize().resolve('conda_locks')
        Files.createDirectories(dir)
        Path out = dir.resolve('batch_correction_env.linux-64.txt')
        if (!Files.exists(out) || out.text != text) {
            Path tmp = Files.createTempFile(dir, 'batch_correction_env', '.tmp')
            tmp.text = text
            Files.move(tmp, out, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE)
        }
        return out.toString()
    }
}
