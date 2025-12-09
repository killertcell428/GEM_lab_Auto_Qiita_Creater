'use client';

import Link from 'next/link';
import StatusBadge from './StatusBadge';
import { format } from 'date-fns';
import { ja } from 'date-fns/locale/ja';
import { ArticleListItem } from '@/lib/api';

interface ArticleStatusCardProps {
  article: ArticleListItem;
}

export default function ArticleStatusCard({ article }: ArticleStatusCardProps) {
  return (
    <Link href={`/articles/${article.id}`}>
      <div className="bg-white dark:bg-gray-800 rounded-lg shadow p-6 hover:shadow-lg transition-shadow cursor-pointer">
        <div className="flex justify-between items-start mb-4">
          <h3 className="text-lg font-semibold text-gray-900 dark:text-white line-clamp-2">
            {article.title}
          </h3>
          <StatusBadge status={article.uiStatusText} phase={article.phase} />
        </div>
        <p className="text-sm text-gray-600 dark:text-gray-400 mb-4">
          {article.nextActionHint}
        </p>
        <div className="flex justify-between items-center text-xs text-gray-500 dark:text-gray-500">
          <span>
            {article.createdAt ? format(new Date(article.createdAt), 'yyyy年MM月dd日 HH:mm', { locale: ja }) : '-'}
          </span>
          <span>
            {article.updatedAt ? format(new Date(article.updatedAt), 'yyyy年MM月dd日 HH:mm', { locale: ja }) : '-'}
          </span>
        </div>
      </div>
    </Link>
  );
}

