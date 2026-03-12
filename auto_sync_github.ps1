$ErrorActionPreference = "Stop"

$repoDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$timestamp = Get-Date -Format "yyyy-MM-dd HH:mm:ss"

Write-Output "[$timestamp] Auto-sync started for $repoDir"

# Pull latest remote changes first to keep local copy current.
git -C $repoDir pull --rebase --autostash origin main
if ($LASTEXITCODE -ne 0) {
    Write-Error "Git pull failed. Auto-sync aborted."
    exit 1
}

# Push any local commits (if branch is ahead).
git -C $repoDir push origin main
if ($LASTEXITCODE -ne 0) {
    Write-Error "Git push failed. Auto-sync aborted."
    exit 1
}

Write-Output "[$timestamp] Auto-sync completed successfully."
