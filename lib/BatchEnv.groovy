import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.security.MessageDigest

/**
 * The batch-correction conda lockfile, ready for Nextflow to build the environment.
 *
 * ConQuR (ivartb/ConQuR_par) and cqrReg 1.2.1 are on no channel the pipeline uses, so they are
 * built from conda_envs/recipes/ and committed under conda_envs/local_channel/. An explicit
 * lockfile needs absolute URLs; the committed lockfile names these two packages under
 * @LOCAL_CHANNEL@, replaced here by file://<projectDir>/conda_envs/local_channel.
 *
 * The checkout path goes into the URL as it is, the same text scripts/setup_batch_env.sh writes,
 * so it may contain only letters, digits and . _ / + - (anything else stops the run).
 * The result is written to <workDir>/conda_locks/ under a name that includes a digest of its
 * content, and never rewritten: runs from different checkouts sharing a work directory get
 * different files, and a file Nextflow has hashed cannot change before the env is built.
 * An existing file with other content (edited or damaged) stops the run instead of being used.
 * The checkout itself is not written here (conda.cacheDir still defaults to it: nextflow.config).
 */
class BatchEnv {
    static final String PLACEHOLDER = '@LOCAL_CHANNEL@'
    static final String LOCKFILE = 'conda_envs/lock/batch_correction_env.linux-64.txt'
    static final String SAFE_PATH = /[A-Za-z0-9._\/+-]+/

    static String lockfile(Object projectDir, Object workDir) {
        Path root = (projectDir as Path).toAbsolutePath().normalize()
        if (!(root.toString() ==~ SAFE_PATH))
            throw new IllegalStateException("The pipeline checkout path may contain only letters, digits and . _ / + - (it is written into the batch lockfile as a file:// URL): ${root}")
        String template = root.resolve(LOCKFILE).text
        if (!template.contains(PLACEHOLDER))
            throw new IllegalStateException("${LOCKFILE} does not reference the local channel")
        String text = template.replace(PLACEHOLDER, 'file://' + root.toString() + '/conda_envs/local_channel')
        String digest = MessageDigest.getInstance('SHA-256').digest(text.getBytes('UTF-8')).encodeHex().toString().substring(0, 16)
        Path dir = (workDir as Path).toAbsolutePath().normalize().resolve('conda_locks')
        Files.createDirectories(dir)
        Path out = dir.resolve('batch_correction_env.linux-64.' + digest + '.txt')
        if (!Files.exists(out)) {
            Path tmp = Files.createTempFile(dir, 'batch_correction_env', '.tmp')
            tmp.setText(text, 'UTF-8')
            Files.move(tmp, out, StandardCopyOption.REPLACE_EXISTING, StandardCopyOption.ATOMIC_MOVE)
        }
        if (out.getText('UTF-8') != text)
            throw new IllegalStateException("${out} does not have the content its name records; it was changed after it was written. Remove it and run again")
        return out.toString()
    }
}
