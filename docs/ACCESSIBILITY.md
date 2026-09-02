# Accessibility and responsive design

G-Synth targets WCAG 2.2 AA for the application-owned interface. Forms use
visible programmatic labels, keyboard-visible focus, semantic landmarks and
live status regions. Status meaning is not conveyed by colour alone, and
application controls provide a minimum 44-pixel target. Password fields have
keyboard-accessible show/hide controls with explicit accessible names.

Scientific sequences, maps and tables use bounded horizontal scrolling when
their information cannot be safely reflowed. The surrounding page remains
responsive, and navigation, authentication, design, cloning and sequencing
validation retain their hierarchy on narrow screens. Reduced-motion mode uses
a static busy indicator while continuing to announce progress.

Automated tests cover core keyboard and semantic behaviour, responsive layout,
password visibility and scientific component rendering. These engineering
checks do not constitute formal WCAG conformance. Release validation should
also include physical-device zoom and assistive-technology testing with the
browser, operating system and observed outcome recorded.
