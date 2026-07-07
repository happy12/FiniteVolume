# `residual_file`/`wall_forces_file` updates appear in bursts, not on the configured interval

## Symptom

A caller polling `residual_file` (or `wall_forces_file`) from outside the
process while a solve is still running -- e.g. Bravo's Airfoil CFD dialog,
which reads it every couple of seconds to redraw the residual chart -- sees
long stretches where the file doesn't change at all, then a jump where many
rows appear at once (e.g. it looked "stuck" until step ~690, then caught up).
This happens regardless of how small `residual_interval`/`wall_forces_interval`
is set to (e.g. 10) -- a small interval does not make updates visible sooner.

## Root cause

`residual_out`/`wall_forces_out` are `std::ofstream`s opened once at the start
of the run and kept open for its entire duration, with one row appended every
`residual_interval`/`wall_forces_interval` steps -- see the identical pattern
repeated in each of the six production run functions:

```cpp
// e.g. src/main.cpp:4172 (run_ransSST) -- same shape in run_diffusion,
// run_advection_diffusion, run_euler, run_navier_stokes, run_ransSA
residual_out.open(case_input.residual_file);
...
if (tracking_residual && step_index % case_input.residual_interval == 0) {
    residual_out << step_index << "," << r.rho << ...;   // src/main.cpp:4246
}
if (tracking_wall_forces && step_index % case_input.wall_forces_interval == 0) {
    write_wall_forces_rows(wall_forces_out, step_index, mesh, reports);   // src/main.cpp:4253
}
```

There is no `.flush()`/`fflush()` anywhere in `main.cpp` (confirmed: zero
matches). Every `<<` above lands in the `ofstream`'s own internal buffer, not
on disk. That buffer is only handed to the OS -- and therefore becomes
visible to any other process reading the same file path -- when it happens to
fill up, or when the stream closes at the run's end. `residual_interval`
controls how often a row is *appended to the buffer*; it has no bearing on
how often that buffer's contents actually *reach the file*. Those are two
independent things, and only the second one determines what an external
poller sees mid-run.

This is distinct from `write_interval`'s numbered VTK/field snapshots and
`wall_profile_interval`'s snapshots (`write_ransSST_fields`,
`write_wall_profile_snapshot`, etc.) -- each of those opens and closes its own
uniquely-numbered file per write, and closing a stream always flushes it, so
those are not affected. Only the two long-lived, run-duration streams
(`residual_out`, `wall_forces_out`) have this gap.

## Recommended fix

Add `residual_out.flush();` immediately after the `residual_out <<` row write,
and `wall_forces_out.flush();` immediately after `write_wall_forces_rows(...)`,
in each of the six write sites (one per run function):

- `run_diffusion` -- around `src/main.cpp:3096`
- `run_advection_diffusion` -- around `src/main.cpp:3219`
- `run_euler` -- around `src/main.cpp:3446`
- `run_navier_stokes` -- around `src/main.cpp:3684`/`3693`
- `run_ransSA` -- around `src/main.cpp:3944`/`3948`
- `run_ransSST` -- around `src/main.cpp:4246`/`4250`

### Why this is worth doing

- **It's the documented purpose of the interval.** `residual_interval`'s own
  doc comment says rows are written "for solution monitoring" (MANUAL.md,
  `CaseInput.h:391`) -- monitoring implies a live external reader, and that
  promise is currently broken: the interval governs buffer writes, not
  visibility.
- **Negligible cost.** `residual_interval`/`wall_forces_interval` already
  throttle how often this code path runs at all (typically every 10+ steps,
  not every step) -- flushing there is nowhere near the pathological
  "flush every single step" case that would actually show up in profiling.
  A flush after an already-infrequent write is cheap.
- **No behavior change for any existing consumer.** The file's final contents
  at run end are identical either way; this only affects how promptly partial
  content becomes visible to a concurrent reader.

## Not in scope for this note

This only covers making the *existing* interval-driven writes visible
promptly. It doesn't propose changing default intervals, adding a new
"flush interval" knob, or touching the numbered-snapshot writers, which
don't need it.
