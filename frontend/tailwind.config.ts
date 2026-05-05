import type { Config } from "tailwindcss";

const config: Config = {
  content: [
    "./app/**/*.{js,ts,jsx,tsx,mdx}",
    "./components/**/*.{js,ts,jsx,tsx,mdx}",
  ],
  theme: {
    extend: {
      colors: {
        canvas: "#f6f7f4",
        ink: "#16211d",
        accent: "#2f6f57",
        muted: "#66756f",
        panel: "#ffffff",
        line: "#d8e0dc",
      },
    },
  },
  plugins: [],
};

export default config;
