# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this directory.

## Purpose

This directory is temporary. It holds systems not yet migrated to `ode_catalog/`, each as an `ODESystemSpec` instance. Add new models to the catalog by following `ode_catalog/ADDING_A_MODEL.md`, not here.

Nothing in `ode_solver/` may hardcode a state name, parameter name, or RHS equation, or import from this directory.
