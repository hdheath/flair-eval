/*
 * Run/source keys shared by input parsing and FLAIR_EVAL. Kept apart from main.nf so it can be
 * included twice without re-including the nf-schema plugin functions.
 */

def safeId(String s) {
    return s.replaceAll(/[^A-Za-z0-9._-]/, '-')
}

// Canonical key of a FLAIR source: runs with equal keys share one FLAIR_ENV task.
def flairSourceKey(Map fl) {
    if (fl.release) {
        return ("release:${fl.release}").toString()
    }
    if (fl.git) {
        return ("git:${fl.git}@${fl.ref}").toString()
    }
    if (fl.env) {
        return ("env:${fl.env.replaceAll('/+$', '')}").toString()
    }
    return ("container:${fl.container}").toString()
}
