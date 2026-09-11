export type IconName =
  | "helix"
  | "plate"
  | "microscope"
  | "target"
  | "scales"
  | "check"
  | "cross"
  | "arrowLeft"
  | "arrowRight"
  | "chevronDown"
  | "search"
  | "book";

type Props = {
  name: IconName;

  size?: number;
  className?: string;

  title?: string;
};


const PATHS: Record<IconName, React.ReactNode> = {

  helix: (
    <>

      <path d="M12 3c6 3 6 6 0 9s-6 6 0 9" />
      <path d="M12 3c-6 3-6 6 0 9s6 6 0 9" />
      <path d="M9 7.2h6M9 16.8h6" opacity="0.6" />
    </>
  ),

  plate: (
    <>
      <circle cx="12" cy="12" r="8.5" />
      <path d="M5.4 6.6A8.5 8.5 0 0 1 17.4 5.6" opacity="0.55" />
      <circle cx="9.6" cy="13.4" r="1.5" />
      <circle cx="14.4" cy="10.6" r="1.1" />
      <circle cx="13.6" cy="15" r="0.9" />
    </>
  ),

  microscope: (
    <>

      <path d="M9.6 2.8h3.2a1 1 0 0 1 1 1V6H8.6V3.8a1 1 0 0 1 1-1Z" />
      <path d="M8.6 6h5.2v4a2.6 2.6 0 0 1-5.2 0Z" />
      <path d="M13.6 8.6A6.6 6.6 0 0 1 14.4 21" />
      <path d="M6.8 17.6h6" opacity="0.6" />
      <path d="M3.4 21h17.2" />
    </>
  ),

  target: (
    <>
      <circle cx="12" cy="12" r="8.4" />
      <circle cx="12" cy="12" r="4.6" opacity="0.55" />
      <circle cx="12" cy="12" r="1.3" fill="currentColor" stroke="none" />
    </>
  ),

  scales: (
    <>
      <path d="M12 4v16M7.5 20h9" />
      <path d="M4 7.6h16" />
      <path d="M4 7.6 1.8 13.2a2.9 2.9 0 0 0 4.4 0Z" />
      <path d="M20 7.6l-2.2 5.6a2.9 2.9 0 0 0 4.4 0Z" />
      <circle cx="12" cy="5.6" r="1.1" fill="currentColor" stroke="none" />
    </>
  ),
  check: <path d="m5 12.6 4.4 4.4L19 7.4" />,
  cross: <path d="M6.4 6.4 17.6 17.6M17.6 6.4 6.4 17.6" />,
  arrowLeft: <path d="M19 12H5m0 0 5.6-5.6M5 12l5.6 5.6" />,
  arrowRight: <path d="M5 12h14m0 0-5.6-5.6M19 12l-5.6 5.6" />,
  chevronDown: <path d="m6.5 9 5.5 5.5L17.5 9" />,
  search: (
    <>
      <circle cx="10.8" cy="10.8" r="6.6" />
      <path d="m15.7 15.7 4.1 4.1" />
    </>
  ),
  book: (
    <>
      <path d="M4 19.5A2.5 2.5 0 0 1 6.5 17H20" />
      <path d="M4 4.5A2.5 2.5 0 0 1 6.5 2H20v20H6.5A2.5 2.5 0 0 1 4 19.5v-15Z" />
      <path d="M8 7h8M8 11h5" opacity="0.55" />
    </>
  ),
};

export default function Icon({ name, size = 20, className, title }: Props) {
  return (
    <svg
      className={`icon ${className ?? ""}`}
      width={size}
      height={size}
      viewBox="0 0 24 24"
      fill="none"
      stroke="currentColor"
      strokeWidth={size < 20 ? 1.9 : 1.6}
      strokeLinecap="round"
      strokeLinejoin="round"
      role={title ? "img" : undefined}
      aria-hidden={title ? undefined : true}
      focusable="false"
    >
      {title && <title>{title}</title>}
      {PATHS[name]}
    </svg>
  );
}
