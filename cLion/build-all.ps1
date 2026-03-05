param(
    [switch]$SkipWin32,
    [switch]$SkipLinux
)

$ErrorActionPreference = "Stop"

$cmakeCmd = Get-Command cmake -ErrorAction SilentlyContinue
if ($cmakeCmd) {
    $cmake = $cmakeCmd.Source
} else {
    $fallbacks = @(
        "$env:LOCALAPPDATA\\Programs\\CLion 2025.3.3\\bin\\cmake\\win\\x64\\bin\\cmake.exe",
        "$env:LOCALAPPDATA\\Programs\\CLion 2025.3\\bin\\cmake\\win\\x64\\bin\\cmake.exe",
        "$env:ProgramFiles\\CMake\\bin\\cmake.exe"
    )

    $cmake = $fallbacks | Where-Object { Test-Path $_ } | Select-Object -First 1
    if (-not $cmake) {
        throw "cmake was not found in PATH or common install locations."
    }
}

function Invoke-CMakeStep {
    param(
        [string]$Label,
        [string[]]$ArgList
    )

    Write-Host "==> $Label"
    & $cmake @ArgList
    if ($LASTEXITCODE -ne 0) {
        throw "Step failed: $Label"
    }
}

if (-not $SkipWin32) {
    Invoke-CMakeStep -Label "Configure release" -ArgList @("--preset", "release")
    Invoke-CMakeStep -Label "Build build-release" -ArgList @("--build", "--preset", "build-release")
}

if ($SkipWin32 -and -not $SkipLinux) {
    Invoke-CMakeStep -Label "Configure release-linux" -ArgList @("--preset", "release-linux")
    Invoke-CMakeStep -Label "Build build-release-linux" -ArgList @("--build", "--preset", "build-release-linux")
} elseif (-not $SkipWin32 -and -not $SkipLinux) {
    Write-Host "Linux release build is triggered by the Release profile build."
}

Write-Host "Done."
